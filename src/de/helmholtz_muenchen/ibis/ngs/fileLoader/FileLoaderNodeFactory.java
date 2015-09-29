package de.helmholtz_muenchen.ibis.ngs.fileLoader;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FileLoader" Node.
 * 
 *
 * @author Tim Jeske
 */
public class FileLoaderNodeFactory 
        extends NodeFactory<FileLoaderNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FileLoaderNodeModel createNodeModel() {
        return new FileLoaderNodeModel();
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
    public NodeView<FileLoaderNodeModel> createNodeView(final int viewIndex,
            final FileLoaderNodeModel nodeModel) {
        return new FileLoaderNodeView(nodeModel);
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
        return new FileLoaderNodeDialog();
    }

}

