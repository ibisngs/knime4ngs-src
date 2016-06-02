package de.helmholtz_muenchen.ibis.ngs.picard;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "PicardTools" Node.
 * 
 *
 * @author 
 */
public class PicardToolsNodeFactory 
        extends NodeFactory<PicardToolsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public PicardToolsNodeModel createNodeModel() {
        return new PicardToolsNodeModel();
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
    public NodeView<PicardToolsNodeModel> createNodeView(final int viewIndex,
            final PicardToolsNodeModel nodeModel) {
        return new HTENodeView<PicardToolsNodeModel>(nodeModel);
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
        return new PicardToolsNodeDialog();
    }

}

