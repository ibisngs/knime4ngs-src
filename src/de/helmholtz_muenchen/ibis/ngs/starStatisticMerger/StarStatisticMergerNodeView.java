package de.helmholtz_muenchen.ibis.ngs.starStatisticMerger;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "StarStatisticMergerNode" Node.
 * 
 *
 * @author Michael Kluge
 */
public class StarStatisticMergerNodeView extends NodeView<StarStatisticMergerNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link StarStatisticMergerNodeModel})
     */
    protected StarStatisticMergerNodeView(final StarStatisticMergerNodeModel nodeModel) {
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

