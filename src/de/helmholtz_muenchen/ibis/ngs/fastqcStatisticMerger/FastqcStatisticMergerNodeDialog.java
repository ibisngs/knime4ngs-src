package de.helmholtz_muenchen.ibis.ngs.fastqcStatisticMerger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeView;

/**
 * Dialog for FastQC Statistic Merger Node 
 * @author Michael Kluge
 *
 */
public class FastqcStatisticMergerNodeDialog extends StatisticMergerNodeView {
	
	@Override
	public String getMergerName() {
		return FastqcStatisticMergerNodeModel.MERGER_NAME;
	}
}

