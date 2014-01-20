package de.helmholtz_muenchen.ibis.ngs.bamloader;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "BAMLoader" Node.
 * 
 *
 * @author 
 */
public class BAMLoaderNodeFactory 
        extends NodeFactory<BAMLoaderNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public BAMLoaderNodeModel createNodeModel() {
        return new BAMLoaderNodeModel();
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
    public NodeView<BAMLoaderNodeModel> createNodeView(final int viewIndex,
            final BAMLoaderNodeModel nodeModel) {
        return new BAMLoaderNodeView(nodeModel);
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
        return new BAMLoaderNodeDialog();
    }

}

