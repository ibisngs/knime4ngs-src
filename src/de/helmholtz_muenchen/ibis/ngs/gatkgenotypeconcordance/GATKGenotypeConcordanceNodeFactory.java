package de.helmholtz_muenchen.ibis.ngs.gatkgenotypeconcordance;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKGenotypeConcordance" Node.
 * 
 *
 * @author Maximilan Hastreiter
 */
public class GATKGenotypeConcordanceNodeFactory 
        extends NodeFactory<GATKGenotypeConcordanceNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKGenotypeConcordanceNodeModel createNodeModel() {
        return new GATKGenotypeConcordanceNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<GATKGenotypeConcordanceNodeModel> createNodeView(final int viewIndex,
            final GATKGenotypeConcordanceNodeModel nodeModel) {
        return new GATKGenotypeConcordanceNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new GATKGenotypeConcordanceNodeDialog();
    }

}

