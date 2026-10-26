package org.pinellolab.dashboard

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpServer
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.net.InetSocketAddress
import java.security.MessageDigest
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.Executors
import java.util.concurrent.ScheduledExecutorService
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicBoolean

/** Pipeline-owned, fail-open lifecycle for the existing advanced dashboard. */
class DashboardObserver implements TraceObserver {
    private static final List<String> FIELDS = [
        'task_id', 'hash', 'process', 'name', 'status', 'exit', 'attempt',
        'submit', 'start', 'complete', 'duration', 'realtime', 'peak_rss',
        'container', 'tag', 'workdir'
    ]

    private Session session
    private Map config = [:]
    private Path outdir
    private Path stateDir
    private Path traceFile
    private Path statusFile
    private Path eventLog
    private Path htmlFile
    private Process monitor
    private HttpServer server
    private ScheduledExecutorService scheduler
    private final AtomicBoolean dirty = new AtomicBoolean(false)
    private final Map<String, Map<String, String>> tasks = new ConcurrentHashMap<>()

    @Override
    void onFlowCreate(Session current) {
        session = current
        // Put plugin options under params to avoid Nextflow's unknown config
        // scope warnings and keep the configuration visible in run provenance.
        Object configured = current.params.get('DASHBOARD_PLUGIN')
        config = configured instanceof Map ? (Map) configured : [:]
        if (!enabled(config.get('enabled'), true)) return
        try {
            outdir = Paths.get(String.valueOf(current.params.outdir)).toAbsolutePath().normalize()
            stateDir = outdir.resolve('pipeline_dashboard/.live')
            traceFile = stateDir.resolve('nextflow_trace.tsv')
            statusFile = stateDir.resolve('status.json')
            eventLog = stateDir.resolve('events.log')
            htmlFile = outdir.resolve('pipeline_dashboard/execution_dashboard.html')
            Files.createDirectories(stateDir)
            putStatus('running', 0)
            flushTrace()
            scheduler = Executors.newSingleThreadScheduledExecutor({ Runnable action ->
                Thread worker = new Thread(action, 'crispr-dashboard-events')
                worker.daemon = true
                return worker
            } as java.util.concurrent.ThreadFactory)
            scheduler.scheduleWithFixedDelay({ flushIfDirty() } as Runnable, 2, 2, TimeUnit.SECONDS)
            startMonitor()
            startServer()
        }
        catch (Exception error) {
            System.err.println("WARN: dashboard plugin initialization failed; pipeline continues: ${error.message}")
        }
    }

    @Override
    void onProcessPending(TaskHandler handler, TraceRecord record) { recordTask(handler, record, 'NEW') }

    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord record) { recordTask(handler, record, 'SUBMITTED') }

    @Override
    void onProcessStart(TaskHandler handler, TraceRecord record) { recordTask(handler, record, 'RUNNING') }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord record) {
        String status = String.valueOf(record?.get('status') ?: '')
        recordTask(handler, record, status.equalsIgnoreCase('FAILED') ? 'FAILED' : 'COMPLETED')
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord record) { recordTask(handler, record, 'CACHED') }

    @Override
    void onFlowError(TaskHandler handler, TraceRecord record) {
        recordTask(handler, record, 'FAILED')
    }

    @Override
    void onFlowComplete() {
        if (stateDir == null) return
        try {
            flushIfDirty()
            putStatus(session?.success ? 'completed' : 'failed', session?.success ? 0 : 1)
            if (monitor != null && !monitor.waitFor(90, TimeUnit.SECONDS)) {
                System.err.println('WARN: dashboard final render exceeded 90 seconds; local snapshot remains available')
                monitor.destroy()
            }
        }
        catch (Exception error) {
            System.err.println("WARN: dashboard finalization failed; pipeline result is unchanged: ${error.message}")
        }
        finally {
            if (scheduler != null) scheduler.shutdownNow()
            if (server != null) server.stop(0)
        }
    }

    private void recordTask(TaskHandler handler, TraceRecord record, String status) {
        if (stateDir == null || handler == null) return
        try {
            String id = String.valueOf(record?.get('task_id') ?: handler.task?.id ?: handler.hashCode())
            Map<String, String> previous = tasks.get(id) ?: [:]
            Map<String, String> row = new LinkedHashMap<>(previous)
            FIELDS.each { String field ->
                Object value = record?.get(field)
                if (value != null && field in ['duration', 'realtime', 'peak_rss']) {
                    value = record.getFmtStr(field)
                }
                if (value != null && String.valueOf(value)) row.put(field, String.valueOf(value))
            }
            row.put('task_id', id)
            row.put('status', status)
            if (!row.get('name')) row.put('name', String.valueOf(handler.task?.name ?: 'unknown'))
            if (!row.get('process')) row.put('process', row.get('name'))
            tasks.put(id, row)
            dirty.set(true)
            Files.writeString(eventLog, "TaskHandler[id: ${id}; name: ${row.get('name')}; status: ${status}]\n",
                StandardCharsets.UTF_8, java.nio.file.StandardOpenOption.CREATE,
                java.nio.file.StandardOpenOption.APPEND)
        }
        catch (Exception error) {
            System.err.println("WARN: dashboard could not record a task event: ${error.message}")
        }
    }

    private void flushIfDirty() {
        if (dirty.getAndSet(false)) flushTrace()
    }

    private synchronized void flushTrace() {
        if (traceFile == null) return
        List<String> lines = [FIELDS.join('\t')]
        tasks.values().sort { a, b -> (a.get('task_id') ?: '') <=> (b.get('task_id') ?: '') }.each { row ->
            lines.add(FIELDS.collect { field -> clean(row.get(field)) }.join('\t'))
        }
        atomicWrite(traceFile, lines.join('\n') + '\n')
    }

    private void putStatus(String status, int exitCode) {
        atomicWrite(statusFile, "{\"status\":\"${status}\",\"exit_code\":${exitCode}}\n")
    }

    private static String clean(String value) {
        return (value ?: '').replaceAll(/[\r\n\t]/, ' ')
    }

    private static void atomicWrite(Path target, String value) {
        Path temp = Files.createTempFile(target.parent, '.dashboard-', '.tmp')
        Files.writeString(temp, value, StandardCharsets.UTF_8)
        try {
            Files.move(temp, target, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
        }
        catch (java.nio.file.AtomicMoveNotSupportedException ignored) {
            Files.move(temp, target, StandardCopyOption.REPLACE_EXISTING)
        }
    }

    private void startMonitor() {
        Object custom = config.get('rendererCommand')
        List<String> command
        if (custom instanceof List && !((List) custom).isEmpty()) {
            command = ((List) custom).collect { String.valueOf(it) }
        }
        else {
            Map resolvedParams = session.config.get('params') instanceof Map ? (Map) session.config.get('params') : [:]
            Map containers = resolvedParams.get('containers') instanceof Map ? (Map) resolvedParams.get('containers') : [:]
            String image = String.valueOf(config.get('image') ?: containers.get('base') ?: '')
            if (!image) throw new IllegalArgumentException('set DASHBOARD_PLUGIN.image or params.containers.base')
            String imageArgument = cachedImageOrUri(image)
            command = ['apptainer', 'exec', '--bind', "${session.baseDir}:${session.baseDir}",
                '--bind', "${outdir}:${outdir}", imageArgument, 'python'].collect { String.valueOf(it) }
        }
        command.addAll([
            session.baseDir.resolve('bin/wandb_html_monitor.py').toString(),
            '--project', String.valueOf(config.get('project') ?: 'crispr-pipeline'),
            '--entity', String.valueOf(config.get('entity') ?: ''),
            '--token-env', String.valueOf(config.get('tokenEnv') ?: 'WB_IGVF'),
            '--run-name', String.valueOf(session.runName ?: session.uniqueId),
            '--source-run-id', String.valueOf(session.uniqueId),
            '--wandb-run-id', String.valueOf(session.uniqueId).replace('-', '').take(20),
            '--replace-run', 'true', '--publish-live-html', 'true',
            '--upload', String.valueOf(enabled(config.get('upload'), true)),
            '--outdir', outdir.toString(), '--trace', traceFile.toString(),
            '--nextflow-log', eventLog.toString(), '--status-file', statusFile.toString(),
            '--dashboard-html', stateDir.resolve('pipeline_execution.html').toString(),
            '--poll-seconds', String.valueOf(config.get('pollSeconds') ?: 20)
        ])
        ProcessBuilder builder = new ProcessBuilder(command)
        builder.redirectErrorStream(true)
        builder.redirectOutput(ProcessBuilder.Redirect.appendTo(stateDir.resolve('monitor.log').toFile()))
        monitor = builder.start()
    }

    private String cachedImageOrUri(String image) {
        String name = image.replaceAll(/[\/:]/, '-') + '.img'
        Object preferred = config.get('imageLibraryDir')
        if (preferred) {
            Path candidate = Paths.get(String.valueOf(preferred)).resolve(name)
            if (Files.isRegularFile(candidate)) return candidate.toString()
        }
        for (String engine : ['apptainer', 'singularity']) {
            Object settings = session.config.get(engine)
            Object library = settings instanceof Map ? ((Map) settings).get('libraryDir') : null
            if (library) {
                Path candidate = Paths.get(String.valueOf(library)).resolve(name)
                if (Files.isRegularFile(candidate)) return candidate.toString()
            }
        }
        return "docker://${image}"
    }

    private void startServer() {
        Map options = config.get('server') instanceof Map ? (Map) config.get('server') : [:]
        if (!enabled(options.get('enabled'), false)) return
        String host = String.valueOf(options.get('host') ?: '127.0.0.1')
        int port = Integer.parseInt(String.valueOf(options.get('port') ?: '8787'))
        if (port < 0 || port > 65535) throw new IllegalArgumentException('dashboard.server.port must be 0–65535')
        boolean local = host in ['127.0.0.1', 'localhost', '::1']
        String token = System.getenv(String.valueOf(options.get('tokenEnv') ?: 'CRISPR_DASHBOARD_TOKEN'))
        if (!local && !token) {
            throw new IllegalArgumentException('non-loopback dashboard binding requires CRISPR_DASHBOARD_TOKEN')
        }
        server = HttpServer.create(new InetSocketAddress(host, port), 0)
        server.createContext('/') { HttpExchange request ->
            if (token && !authorized(request, token)) {
                request.responseHeaders.add('WWW-Authenticate', 'Basic realm="CRISPR dashboard"')
                respond(request, 401, 'Authentication required'.getBytes(StandardCharsets.UTF_8), 'text/plain')
                return
            }
            if (request.requestURI.path == '/health') {
                respond(request, 200, 'ok\n'.getBytes(StandardCharsets.UTF_8), 'text/plain')
                return
            }
            if (!(request.requestURI.path in ['/', '/dashboard'])) {
                respond(request, 404, 'Not found'.getBytes(StandardCharsets.UTF_8), 'text/plain')
                return
            }
            if (!Files.isRegularFile(htmlFile)) {
                respond(request, 503, 'Dashboard is starting; refresh shortly.'.getBytes(StandardCharsets.UTF_8), 'text/plain')
                return
            }
            String html = Files.readString(htmlFile, StandardCharsets.UTF_8)
            int refresh = Integer.parseInt(String.valueOf(options.get('refreshSeconds') ?: '30'))
            if (refresh > 0) html = html.replaceFirst('</head>', "<meta http-equiv=\"refresh\" content=\"${refresh}\"></head>")
            respond(request, 200, html.getBytes(StandardCharsets.UTF_8), 'text/html; charset=utf-8')
        }
        server.executor = Executors.newCachedThreadPool({ Runnable action ->
            Thread worker = new Thread(action, 'crispr-dashboard-http')
            worker.daemon = true
            return worker
        } as java.util.concurrent.ThreadFactory)
        server.start()
        System.err.println("CRISPR dashboard: http://${host}:${server.address.port}/ (available while Nextflow runs)")
    }

    private static boolean authorized(HttpExchange request, String token) {
        String header = request.requestHeaders.getFirst('Authorization') ?: ''
        if (!header.startsWith('Basic ')) return false
        try {
            String credentials = new String(Base64.decoder.decode(header.substring(6)), StandardCharsets.UTF_8)
            byte[] observed = credentials.getBytes(StandardCharsets.UTF_8)
            byte[] expected = "dashboard:${token}".getBytes(StandardCharsets.UTF_8)
            return MessageDigest.isEqual(observed, expected)
        }
        catch (Exception ignored) { return false }
    }

    private static void respond(HttpExchange request, int status, byte[] body, String type) {
        request.responseHeaders.add('Content-Type', type)
        request.responseHeaders.add('Cache-Control', 'no-store')
        request.responseHeaders.add('X-Content-Type-Options', 'nosniff')
        request.sendResponseHeaders(status, body.length)
        request.responseBody.withCloseable { it.write(body) }
    }

    private static boolean enabled(Object value, boolean fallback) {
        if (value == null) return fallback
        return String.valueOf(value).toLowerCase() in ['true', '1', 'yes']
    }
}
