package de.helmholtz_muenchen.ibis.ngs.featureCountsStatisticMerger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeDialog;

/**
 * Dialog for FeatureCount Statistic Merger Node 
 * @author Michael Kluge
 *
 */
public class FeatureCountsStatisticMergerNodeDialog extends StatisticMergerNodeDialog {
	
	@Override
	public String getMergerName() {
		return FeatureCountsStatisticMergerNodeModel.MERGER_NAME;
	}
}