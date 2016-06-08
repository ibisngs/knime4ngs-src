package de.helmholtz_muenchen.ibis.ngs.VCFSorter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "VCFSorter" Node.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class VCFSorterNodeFactory 
        extends NodeFactory<VCFSorterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VCFSorterNodeModel createNodeModel() {
        return new VCFSorterNodeModel();
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
    public NodeView<VCFSorterNodeModel> createNodeView(final int viewIndex,
            final VCFSorterNodeModel nodeModel) {
        return new HTENodeView<VCFSorterNodeModel>(nodeModel);
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
        return new VCFSorterNodeDialog();
    }

}

