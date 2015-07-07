package de.helmholtz_muenchen.ibis.ngs.samtoolshybrid;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SamtoolsHybrid" Node.
 * 
 *
 * @author Tanzeem Haque
 */
public class SamtoolsHybridNodeFactory 
        extends NodeFactory<SamtoolsHybridNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SamtoolsHybridNodeModel createNodeModel() {
        return new SamtoolsHybridNodeModel();
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
    public NodeView<SamtoolsHybridNodeModel> createNodeView(final int viewIndex,
            final SamtoolsHybridNodeModel nodeModel) {
        return new SamtoolsHybridNodeView(nodeModel);
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
        return new SamtoolsHybridNodeDialog();
    }

}

