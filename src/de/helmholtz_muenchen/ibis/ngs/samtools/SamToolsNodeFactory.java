package de.helmholtz_muenchen.ibis.ngs.samtools;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SamTools" Node.
 * 
 *
 * @author 
 */
public class SamToolsNodeFactory 
        extends NodeFactory<SamToolsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SamToolsNodeModel createNodeModel() {
        return new SamToolsNodeModel();
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
    public NodeView<SamToolsNodeModel> createNodeView(final int viewIndex,
            final SamToolsNodeModel nodeModel) {
        return new SamToolsNodeView(nodeModel);
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
        return new SamToolsNodeDialog();
    }

}

