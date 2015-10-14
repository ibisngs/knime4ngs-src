package de.helmholtz_muenchen.ibis.ngs.geneticBackgroundModel;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GeneticBackgroundModel" Node.
 * 
 *
 * @author Tim Jeske
 */
public class GeneticBackgroundModelNodeFactory 
        extends NodeFactory<GeneticBackgroundModelNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GeneticBackgroundModelNodeModel createNodeModel() {
        return new GeneticBackgroundModelNodeModel();
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
    public NodeView<GeneticBackgroundModelNodeModel> createNodeView(final int viewIndex,
            final GeneticBackgroundModelNodeModel nodeModel) {
        return new GeneticBackgroundModelNodeView(nodeModel);
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
        return new GeneticBackgroundModelNodeDialog();
    }

}

