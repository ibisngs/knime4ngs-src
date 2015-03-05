package de.helmholtz_muenchen.ibis.ngs.gatkcombinegvcfs;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "CombineGVCFs" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class CombineGVCFsNodeFactory 
        extends NodeFactory<CombineGVCFsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public CombineGVCFsNodeModel createNodeModel() {
        return new CombineGVCFsNodeModel();
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
    public NodeView<CombineGVCFsNodeModel> createNodeView(final int viewIndex,
            final CombineGVCFsNodeModel nodeModel) {
        return new CombineGVCFsNodeView(nodeModel);
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
        return new CombineGVCFsNodeDialog();
    }

}
