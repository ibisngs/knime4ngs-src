package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulator;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "RawReadManipulator" Node.
 * 
 *
 * @author 
 */
public class RawReadManipulatorNodeFactory 
        extends NodeFactory<RawReadManipulatorNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public RawReadManipulatorNodeModel createNodeModel() {
        return new RawReadManipulatorNodeModel();
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
    public NodeView<RawReadManipulatorNodeModel> createNodeView(final int viewIndex,
            final RawReadManipulatorNodeModel nodeModel) {
        return new HTENodeView<RawReadManipulatorNodeModel>(nodeModel);
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
        return new RawReadManipulatorNodeDialog();
    }

}

