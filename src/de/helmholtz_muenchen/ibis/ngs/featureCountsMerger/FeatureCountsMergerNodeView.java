package de.helmholtz_muenchen.ibis.ngs.featureCountsMerger;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FeatureCountsMerger" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FeatureCountsMergerNodeView extends NodeView<FeatureCountsMergerNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FeatureCountsMergerNodeModel})
     */
    protected FeatureCountsMergerNodeView(final FeatureCountsMergerNodeModel nodeModel) {
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

