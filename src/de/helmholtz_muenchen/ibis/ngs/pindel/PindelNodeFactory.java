package de.helmholtz_muenchen.ibis.ngs.pindel;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Pindel" Node.
 * 
 *
 * @author 
 */
public class PindelNodeFactory 
        extends NodeFactory<PindelNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public PindelNodeModel createNodeModel() {
        return new PindelNodeModel();
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
    public NodeView<PindelNodeModel> createNodeView(final int viewIndex,
            final PindelNodeModel nodeModel) {
        return new PindelNodeView(nodeModel);
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
        return new PindelNodeDialog();
    }

}

