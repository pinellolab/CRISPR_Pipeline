package org.pinellolab.dashboard

import nextflow.Session
import spock.lang.Specification

class DashboardObserverFactoryTest extends Specification {
    def 'registers one dashboard trace observer'() {
        when:
        def observers = new DashboardObserverFactory().create(Mock(Session))

        then:
        observers.size() == 1
        observers.first() instanceof DashboardObserver
    }
}
