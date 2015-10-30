package de.helmholtz_muenchen.ibis.ngs.featureCounts;

import org.knime.core.node.NodeView;


/**
 * <code>NodeView</code> for the "FeatureCounts" Node.
 * featureCounts: an efficient general purpose program for assigning sequence reads to genomic features
 * For details, please see paper: 
 * Liao et al, Bioinformatics 2013; doi: 10.1093/bioinformatics/btt656
 *
 * @author Michael Kluge
 */
public class FeatureCountsNodeView extends NodeView<FeatureCountsNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FeatureCountsNodeModel})
     */
    protected FeatureCountsNodeView(final FeatureCountsNodeModel nodeModel) {
        super(nodeModel);
    }

	@Override
	protected void onClose() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void onOpen() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void modelChanged() {
		// TODO Auto-generated method stub
		
	}
}

