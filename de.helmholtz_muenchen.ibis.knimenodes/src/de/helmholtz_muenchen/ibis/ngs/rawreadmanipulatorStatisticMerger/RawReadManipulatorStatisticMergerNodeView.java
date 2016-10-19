package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulatorStatisticMerger;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "RawReadManipulatorStatisticMergerNode" Node.
 * 
 *
 * @author Michael Kluge
 */
public class RawReadManipulatorStatisticMergerNodeView extends NodeView<RawReadManipulatorStatisticMergerNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link RawReadManipulatorStatisticMergerNodeModel})
     */
    protected RawReadManipulatorStatisticMergerNodeView(final RawReadManipulatorStatisticMergerNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
    }
}

