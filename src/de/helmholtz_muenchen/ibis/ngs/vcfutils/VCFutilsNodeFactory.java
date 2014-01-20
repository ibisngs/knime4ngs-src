package de.helmholtz_muenchen.ibis.ngs.vcfutils;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "VCFutils" Node.
 * 
 *
 * @author 
 */
public class VCFutilsNodeFactory 
        extends NodeFactory<VCFutilsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VCFutilsNodeModel createNodeModel() {
        return new VCFutilsNodeModel();
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
    public NodeView<VCFutilsNodeModel> createNodeView(final int viewIndex,
            final VCFutilsNodeModel nodeModel) {
        return new VCFutilsNodeView(nodeModel);
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
        return new VCFutilsNodeDialog();
    }

}

