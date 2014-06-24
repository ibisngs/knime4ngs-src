package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "MatsResultIndexer" Node.
 * 
 *
 * @author Michael Kluge
 */
public class MatsResultIndexerNodeFactory 
        extends NodeFactory<MatsResultIndexerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public MatsResultIndexerNodeModel createNodeModel() {
        return new MatsResultIndexerNodeModel();
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
    public NodeView<MatsResultIndexerNodeModel> createNodeView(final int viewIndex,
            final MatsResultIndexerNodeModel nodeModel) {
        return new MatsResultIndexerNodeView(nodeModel);
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
        return new MatsResultIndexerNodeDialog();
    }

}

