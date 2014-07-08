package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulatorStatisticMerger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeDialog;

/**
 * Dialog for RawReadManipulator Statistic Merger Node 
 * @author Michael Kluge
 *
 */
public class RawReadManipulatorStatisticMergerNodeDialog extends StatisticMergerNodeDialog {
	
	@Override
	public String getMergerName() {
		return RawReadManipulatorStatisticMergerNodeModel.MERGER_NAME;
	}
}