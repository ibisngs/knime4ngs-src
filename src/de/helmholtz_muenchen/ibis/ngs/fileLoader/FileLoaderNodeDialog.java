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
package de.helmholtz_muenchen.ibis.ngs.fileLoader;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;


/**
 * <code>NodeDialog</code> for the "FileLoader" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class FileLoaderNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString m_infile1 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE1,"");
	private final SettingsModelString m_infile2 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE2,"");
	private final SettingsModelBoolean m_isList = new SettingsModelBoolean(FileLoaderNodeModel.CFGKEY_ISLIST,false);

	
	/**
     * New pane for configuring the FileLoader node.
     */
    protected FileLoaderNodeDialog() {
    	addDialogComponent(new DialogComponentBoolean(m_isList,"Handle input file as list of file paths?"));
    	
    	createNewGroup("Input file (BAM/SAM, fastQ, VCF)");
    	addDialogComponent(new DialogComponentFileChooser(m_infile1, "his0_infile1", 0, ".fa|.fasta",".bam|.sam","fastq|.fastq.gz|.fq|.fq.gz",".vcf|.vcf.gz|.gvcf",".csv|.tsv|.list"));

    	createNewGroup("Second optional fastQ file");
    	addDialogComponent(new DialogComponentFileChooser(m_infile2, "his0_infile2", 0, "fastq|.fastq.gz|.fq"));
    
    	
    	m_infile1.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				String file1 = m_infile1.getStringValue();
				file1 = IO.removeZipExtension(file1);
				if((file1.endsWith(".fastq") || file1.endsWith(".fq")) && !m_isList.getBooleanValue()) {
					m_infile2.setEnabled(true);
				} else {
					m_infile2.setEnabled(false);
					m_infile2.setStringValue("");
				}
			}
    	});
    	
    	m_isList.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(m_isList.getBooleanValue()) {
					m_infile2.setEnabled(false);
				} else {
					String file1 = m_infile1.getStringValue();
					file1 = IO.removeZipExtension(file1);
					if((file1.endsWith(".fastq") || file1.endsWith(".fq")) && !m_isList.getBooleanValue()) {
						m_infile2.setEnabled(true);
					}
				}
			}
    	});
    
    }
    
}

