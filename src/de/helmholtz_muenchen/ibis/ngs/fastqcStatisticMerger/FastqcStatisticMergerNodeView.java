package de.helmholtz_muenchen.ibis.ngs.fastqcStatisticMerger;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FastqcStatisticMergerNodeModel" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FastqcStatisticMergerNodeView extends NodeView<FastqcStatisticMergerNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FastqcStatisticMergerNodeModel})
     */
    protected FastqcStatisticMergerNodeView(final FastqcStatisticMergerNodeModel nodeModel) {
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

