package de.helmholtz_muenchen.ibis.ngs.gatkCombineVCFs;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "CombineVCFs" Node.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class CombineVCFsNodeFactory 
        extends NodeFactory<CombineVCFsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public CombineVCFsNodeModel createNodeModel() {
        return new CombineVCFsNodeModel();
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
    public NodeView<CombineVCFsNodeModel> createNodeView(final int viewIndex,
            final CombineVCFsNodeModel nodeModel) {
        return new HTENodeView<CombineVCFsNodeModel>(nodeModel);
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
        return new CombineVCFsNodeDialog();
    }

}

