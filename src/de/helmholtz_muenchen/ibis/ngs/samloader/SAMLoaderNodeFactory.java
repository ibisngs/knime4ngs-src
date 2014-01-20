package de.helmholtz_muenchen.ibis.ngs.samloader;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SAMLoader" Node.
 * 
 *
 * @author 
 */
public class SAMLoaderNodeFactory 
        extends NodeFactory<SAMLoaderNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SAMLoaderNodeModel createNodeModel() {
        return new SAMLoaderNodeModel();
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
    public NodeView<SAMLoaderNodeModel> createNodeView(final int viewIndex,
            final SAMLoaderNodeModel nodeModel) {
        return new SAMLoaderNodeView(nodeModel);
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
        return new SAMLoaderNodeDialog();
    }

}

