package de.helmholtz_muenchen.ibis.ngs.featureCountsStatisticMerger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeView;

/**
 * Dialog for FeatureCount Statistic Merger Node 
 * @author Michael Kluge
 *
 */
public class FeatureCountsStatisticMergerNodeDialog extends StatisticMergerNodeView {
	
	@Override
	public String getMergerName() {
		return FeatureCountsStatisticMergerNodeModel.MERGER_NAME;
	}
}