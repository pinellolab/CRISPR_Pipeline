package org.pinellolab.dashboard

import nextflow.trace.TraceRecord
import spock.lang.Specification
import spock.lang.TempDir
import java.nio.file.Path

class DashboardObserverTaskTest extends Specification {
    @TempDir Path directory

    def 'an existing absolute image is not rewritten as a Docker URI'() {
        given:
        def image = directory.resolve('pinned.img')
        image.text = 'fixture'
        def method = DashboardObserver.getDeclaredMethod('cachedImageOrUri', String)
        method.accessible = true

        expect:
        method.invoke(new DashboardObserver(), image.toString()) == image.toString()
        method.invoke(new DashboardObserver(), 'docker://ghcr.io/example/image:tag') == 'docker://ghcr.io/example/image:tag'
    }

    def 'cached records without handlers and separate retry attempts keep work paths'() {
        given:
        def observer = new DashboardObserver()
        ['stateDir': directory, 'traceFile': directory.resolve('trace.tsv'),
         'eventLog': directory.resolve('events.log')].each { name, value ->
            def field = DashboardObserver.getDeclaredField(name)
            field.accessible = true
            field.set(observer, value)
        }
        def failed = new TraceRecord()
        failed.put('task_id', 1)
        failed.put('name', 'mapper (sample)')
        failed.put('attempt', 1)
        failed.put('workdir', '/work/old')
        failed.put('status', 'FAILED')
        def cached = new TraceRecord()
        cached.put('task_id', 1)
        cached.put('name', 'mapper (sample)')
        cached.put('attempt', 2)
        cached.put('workdir', '/work/new')

        when:
        observer.onProcessComplete(null, failed)
        observer.onProcessCached(null, cached)
        def flush = DashboardObserver.getDeclaredMethod('flushTrace')
        flush.accessible = true
        flush.invoke(observer)
        def text = directory.resolve('trace.tsv').text

        then:
        text.contains('/work/old')
        text.contains('/work/new')
        text.contains('FAILED')
        text.contains('CACHED')
        text.readLines().size() == 3
    }
}
