package de.helmholtz_muenchen.ibis.ngs.starStatisticMerger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeView;

/**
 * Dialog for Star Statistic Merger Node 
 * @author Michael Kluge
 *
 */
public class StarStatisticMergerNodeDialog extends StatisticMergerNodeView {
	
	@Override
	public String getMergerName() {
		return StarStatisticMergerNodeModel.MERGER_NAME;
	}
}