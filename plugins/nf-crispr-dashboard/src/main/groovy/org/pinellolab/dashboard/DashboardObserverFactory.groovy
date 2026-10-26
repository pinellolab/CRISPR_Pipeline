package org.pinellolab.dashboard

import nextflow.Session
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

class DashboardObserverFactory implements TraceObserverFactory {
    @Override
    Collection<TraceObserver> create(Session session) {
        return [new DashboardObserver()]
    }
}
