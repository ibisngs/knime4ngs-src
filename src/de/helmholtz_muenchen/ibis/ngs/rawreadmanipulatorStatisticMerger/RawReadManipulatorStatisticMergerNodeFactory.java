package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulatorStatisticMerger;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "RawReadManipulatorStatisticMergerNode" Node.
 * 
 *
 * @author Michael Kluge
 */
public class RawReadManipulatorStatisticMergerNodeFactory 
        extends NodeFactory<RawReadManipulatorStatisticMergerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public RawReadManipulatorStatisticMergerNodeModel createNodeModel() {
        return new RawReadManipulatorStatisticMergerNodeModel();
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
    public NodeView<RawReadManipulatorStatisticMergerNodeModel> createNodeView(final int viewIndex, final RawReadManipulatorStatisticMergerNodeModel nodeModel) {
        return new RawReadManipulatorStatisticMergerNodeView(nodeModel);
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
        return new RawReadManipulatorStatisticMergerNodeDialog();
    }

}

