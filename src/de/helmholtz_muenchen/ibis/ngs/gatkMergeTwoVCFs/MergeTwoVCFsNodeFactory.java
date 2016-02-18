package de.helmholtz_muenchen.ibis.ngs.gatkMergeTwoVCFs;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "MergeTwoVCFs" Node.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class MergeTwoVCFsNodeFactory 
        extends NodeFactory<MergeTwoVCFsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public MergeTwoVCFsNodeModel createNodeModel() {
        return new MergeTwoVCFsNodeModel();
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
    public NodeView<MergeTwoVCFsNodeModel> createNodeView(final int viewIndex,
            final MergeTwoVCFsNodeModel nodeModel) {
        return new MergeTwoVCFsNodeView(nodeModel);
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
        return new MergeTwoVCFsNodeDialog();
    }

}

