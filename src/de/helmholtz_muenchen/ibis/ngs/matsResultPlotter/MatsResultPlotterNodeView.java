package de.helmholtz_muenchen.ibis.ngs.matsResultPlotter;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "MatsResultPlotter" Node.
 * 
 *
 * @author Michael Kluge
 */
public class MatsResultPlotterNodeView extends NodeView<MatsResultPlotterNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link MatsResultPlotterNodeModel})
     */
    protected MatsResultPlotterNodeView(final MatsResultPlotterNodeModel nodeModel) {
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

