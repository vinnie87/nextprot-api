package org.nextprot.api.core.utils.graph;

public class CvTermGraphTest extends AbstractCvTermGraphTest {

    @Override
    protected DirectedGraph createGraph() {

        return new CvTermGraph();
    }
}