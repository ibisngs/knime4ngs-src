package de.helmholtz_muenchen.ibis.ngs.bcftools;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Bcftools" Node.
 * 
 *
 * @author Max
 */
public class BcftoolsNodeFactory 
        extends NodeFactory<BcftoolsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public BcftoolsNodeModel createNodeModel() {
        return new BcftoolsNodeModel();
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
    public NodeView<BcftoolsNodeModel> createNodeView(final int viewIndex,
            final BcftoolsNodeModel nodeModel) {
        return new BcftoolsNodeView(nodeModel);
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
        return new BcftoolsNodeDialog();
    }

}

