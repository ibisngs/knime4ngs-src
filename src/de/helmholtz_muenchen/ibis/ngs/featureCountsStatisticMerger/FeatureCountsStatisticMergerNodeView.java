package de.helmholtz_muenchen.ibis.ngs.featureCountsStatisticMerger;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FeatureCountsStatisticMergerNode" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FeatureCountsStatisticMergerNodeView extends NodeView<FeatureCountsStatisticMergerNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FeatureCountsStatisticMergerNodeModel})
     */
    protected FeatureCountsStatisticMergerNodeView(final FeatureCountsStatisticMergerNodeModel nodeModel) {
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

