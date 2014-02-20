package de.helmholtz_muenchen.ibis.ngs.samSelector;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector.FileSelectorNodeModel;


/**
 * <code>FastaSelectorNodeModel</code> can be used to select multiple fasta files.
 *
 * @author Michael Kluge
 */
public class SamSelectorNodeModel extends FileSelectorNodeModel {

	public static final String OUTPUT_NAME_SAM_FILES 	= "SamFiles";
	public static final String FILTERTYPE_NAME 			= "SAM/BAM";
	public static final String FILENAME_END_REGEX 		= "\\.(sam|bam)$";
	
	@Override
	public String getNameOfOutputCol() {
		return OUTPUT_NAME_SAM_FILES;
	}
	
	@Override
	protected String getFiletypeName() {
		return FILTERTYPE_NAME;
	}
}

