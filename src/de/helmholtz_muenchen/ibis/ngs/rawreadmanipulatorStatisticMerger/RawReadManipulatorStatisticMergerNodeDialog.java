package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulatorStatisticMerger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeView;

/**
 * Dialog for RawReadManipulator Statistic Merger Node 
 * @author Michael Kluge
 *
 */
public class RawReadManipulatorStatisticMergerNodeDialog extends StatisticMergerNodeView {
	
	@Override
	public String getMergerName() {
		return RawReadManipulatorStatisticMergerNodeModel.MERGER_NAME;
	}
}