package org.nextprot.api.core.utils.graph;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.google.common.base.Preconditions;
import org.nextprot.api.commons.constants.TerminologyCv;
import org.nextprot.api.commons.utils.graph.DirectedGraph;
import org.nextprot.api.core.domain.CvTerm;
import org.nextprot.api.core.service.TerminologyService;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;
import java.util.logging.Logger;

/**
 * Backing graph of cv terms as ints node as the maximum id is 225719 !!
 *
 * Created by fnikitin on 23.06.17.
 */
public abstract class BaseCvTermGraph implements Serializable {

    private static final long serialVersionUID = 1L;

    private final static Logger LOGGER = Logger.getLogger(BaseCvTermGraph.class.getSimpleName());

    protected final DirectedGraph graph;

    public BaseCvTermGraph(TerminologyCv terminologyCv, TerminologyService service, Supplier<DirectedGraph> graphSupplier) {

        Preconditions.checkNotNull(terminologyCv);
        Preconditions.checkNotNull(service);
        Preconditions.checkNotNull(graphSupplier);

        List<CvTerm> cvTerms = service.findCvTermsByOntology(terminologyCv.name());

        graph = graphSupplier.get();
        graph.setGraphLabel(terminologyCv.name());

        cvTerms.forEach(this::addCvTermNode);
        cvTerms.forEach(this::addCvTermEdges);
    }

    BaseCvTermGraph(TerminologyCv terminologyCv, DirectedGraph graph) {

        Preconditions.checkNotNull(graph);

        this.graph = graph;
        graph.setGraphLabel(terminologyCv.name());
    }

    private void addCvTermNode(CvTerm cvTerm) {

        int nodeId = Math.toIntExact(cvTerm.getId());

        graph.addNode(nodeId);
        graph.setNodeLabel(nodeId, cvTerm.getAccession());
    }

    private void addCvTermEdges(CvTerm cvTerm) {

        List<String> parentAccessions = cvTerm.getAncestorAccession();

        if (parentAccessions != null) {
            parentAccessions.forEach(parent -> {
                try {
                    graph.addEdge(getCvTermIdByAccession(parent), Math.toIntExact(cvTerm.getId()));
                } catch (NotFoundNodeException e) {
                    LOGGER.warning(cvTerm.getAccession()+" cannot connect to unknown node parent: "+e.getMessage());
                }
            });
        }
    }

    /**
     * @return the TerminologyCv of this graph of CvTerms
     */
    public TerminologyCv getTerminologyCv() {
        return TerminologyCv.getTerminologyOf(graph.getGraphLabel());
    }

    /**
     * @return the CvTerm with given id
     */
    public String getCvTermAccessionById(int id) {

        return graph.getNodeLabel(id);
    }

    /**
     * @return the id of cvterm with given accession
     */
    public int getCvTermIdByAccession(String accession) throws NotFoundNodeException {

        int cvTerm = graph.getNode(accession);

        if (cvTerm == -1)
            throw this.new NotFoundNodeException(accession);

        return cvTerm;
    }

    /**
     * @return true if cvTermAccession was found
     */
    public boolean hasCvTermAccession(String cvTermAccession) {

        return graph.getNode(cvTermAccession) != -1;
    }

    int[] getSources() {

        return graph.getSources();
    }

    public int[] getNodes() {

        return graph.getNodes();
    }

    public int countNodes() {

        return graph.countNodes();
    }

    public int countEdges() {

        return graph.countEdges();
    }

    public boolean isAncestorOf(int queryAncestor, int queryDescendant) {

        return graph.isAncestorOf(queryAncestor, queryDescendant);
    }

    public boolean isDescendantOf(int queryDescendant, int queryAncestor) {

        return graph.isDescendantOf(queryDescendant, queryAncestor);
    }

    public int[] getAncestors(int cvTermId) {

        return graph.getAncestors(cvTermId);
    }

    public BaseCvTermGraph calcAncestorSubgraph(int cvTermId) {

        return newSupplier(getTerminologyCv(), graph.calcAncestorSubgraph(cvTermId)).get();
    }

    protected abstract Supplier<? extends BaseCvTermGraph> newSupplier(TerminologyCv terminologyCv, DirectedGraph graph);

    public int[] getParents(int cvTermId) {

        return graph.getPredecessors(cvTermId);
    }

    public int[] getChildren(int cvTermId) {

        return graph.getSuccessors(cvTermId);
    }

    public int[] getEdges() {
        return graph.getEdges();
    }

    public int getTailNode(int edge) {
        return graph.getTailNode(edge);
    }

    public int getHeadNode(int edge) {
        return graph.getHeadNode(edge);
    }

    public View toView() {

        View view = new View();

        view.setLabel(graph.getGraphLabel());

        for (int nid : graph.getNodes()) {
            View.Node node = new View.Node();
            node.setData(nid, graph.getNodeLabel(nid));
            view.addNode(node);
        }

        for (int eid : graph.getEdges()) {
            View.Edge edge = new View.Edge();
            edge.setData(eid, graph.getTailNode(eid), graph.getHeadNode(eid));
            edge.setLabel(graph.getEdgeLabel(eid));
            view.addEdge(edge);
        }

        return view;
    }

    /**
     * Thrown if no nodes map the given cvterm accession
     */
    public class NotFoundNodeException extends Exception {

        public NotFoundNodeException(String accession) {

            super("CvTerm node with accession "+accession+" was not found in "+getTerminologyCv() + " graph");
        }
    }

    public static class View {

        private String label;
        private List<Node> nodes = new ArrayList<>();
        private List<Edge> edges = new ArrayList<>();

        public String getLabel() {
            return label;
        }

        public void setLabel(String label) {
            this.label = label;
        }

        public List<Node> getNodes() {
            return nodes;
        }

        public void addNode(Node node) {
            this.nodes.add(node);
        }

        public List<Edge> getEdges() {
            return edges;
        }

        public void addEdge(Edge edge) {
            this.edges.add(edge);
        }

        public static class Node {

            private int id;
            private String label;

            public int getId() {
                return id;
            }

            public String getLabel() {
                return label;
            }

            public void setData(int id, String label) {
                this.id = id;
                this.label = label;
            }
        }

        @JsonInclude(JsonInclude.Include.NON_NULL)
        public static class Edge {

            private int id;
            private int tail;
            private int head;
            private String label;

            public int getId() {
                return id;
            }

            public void setId(int id) {
                this.id = id;
            }

            public int getTail() {
                return tail;
            }

            public int getHead() {
                return head;
            }

            public void setData(int id, int tail, int head) {
                this.id = id;
                this.tail = tail;
                this.head = head;
            }

            public String getLabel() {
                return label;
            }

            public void setLabel(String label) {
                this.label = label;
            }
        }

    }
}
