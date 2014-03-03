package de.helmholtz_muenchen.ibis.ngs.vcfloader;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "VCFLoader" Node.
 * 
 *
 * @author 
 */
public class VCFLoaderNodeFactory 
        extends NodeFactory<VCFLoaderNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VCFLoaderNodeModel createNodeModel() {
        return new VCFLoaderNodeModel();
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
    public NodeView<VCFLoaderNodeModel> createNodeView(final int viewIndex,
            final VCFLoaderNodeModel nodeModel) {
        return new VCFLoaderNodeView(nodeModel);
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
        return new VCFLoaderNodeDialog();
    }

}

