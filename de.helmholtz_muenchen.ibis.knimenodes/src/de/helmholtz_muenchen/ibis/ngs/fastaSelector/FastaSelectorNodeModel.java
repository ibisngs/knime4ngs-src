/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package de.helmholtz_muenchen.ibis.ngs.fastaSelector;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector.FileSelectorNodeModel;


/**
 * <code>FastaSelectorNodeModel</code> can be used to select multiple fasta files.
 *
 * @author Michael Kluge
 */
public class FastaSelectorNodeModel extends FileSelectorNodeModel {

	public static final String OUTPUT_NAME_FASTA_FILES 	= "FastaFiles";
	public static final String FILTERTYPE_NAME 			= "FASTA";
	public static final String FILENAME_END_REGEX 		= "\\.(fa|fasta)$";
	
	@Override
	public String getNameOfOutputCol() {
		return OUTPUT_NAME_FASTA_FILES;
	}
	
	@Override
	protected String getFiletypeName() {
		return FILTERTYPE_NAME;
	}
}

