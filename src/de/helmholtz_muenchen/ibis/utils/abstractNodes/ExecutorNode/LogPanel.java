package de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridLayout;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;

public class LogPanel extends JTabbedPane {
	private static final long serialVersionUID = 8322929670792162817L;

	String STDOUT;
	String STDERR;

	JPanel panel_stdout;
	JPanel panel_stderr;
	
	JTextArea textarea_stdout;
	JTextArea textarea_stderr;

	JScrollPane scollpane_stdout;
	JScrollPane scollpane_stderr;
	
	public LogPanel(String stdout, String stderr){
		super();

		STDOUT = stdout;
		STDERR = stderr;

		setPreferredSize(new Dimension(800, 600));

		init();
	}
	
	private void init(){
		// Textareas
		textarea_stdout = new JTextArea();
		textarea_stdout.setLineWrap(false);
		textarea_stdout.setEditable(false);
		
		textarea_stderr = new JTextArea();
		textarea_stderr.setLineWrap(false);
		textarea_stderr.setEditable(false);
		
		// scrollpanes
		scollpane_stdout = new JScrollPane(textarea_stdout, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		scollpane_stderr = new JScrollPane(textarea_stderr, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		
		// panels
		panel_stdout = new JPanel(new GridLayout(1,1));
		panel_stdout.add(scollpane_stdout);
		
		panel_stderr = new JPanel(new GridLayout(1,1));
		panel_stderr.add(scollpane_stderr);

		// tabs
		this.addTab("STDOUT", panel_stdout);
		this.addTab("STDERR", panel_stderr);

		// add text and update
		this.updateView();
	}
	
	private void setText(){
		if(STDOUT == null){
			textarea_stdout.setText(ExecutorNodeModel.LOGMESSAGE_LOG_DISABLED);	
		}else{
			textarea_stdout.setText(STDOUT);	
		}
		if(STDERR == null){
			textarea_stderr.setText(ExecutorNodeModel.LOGMESSAGE_LOG_DISABLED);	
		}else{
			textarea_stderr.setText(STDERR);	
		}
	}
	/**
	 * @see javax.swing.JComponent#paint(java.awt.Graphics)
	 */
	@Override
	public void paint(Graphics g) {
		super.paint(g);
		setText();
	}

	/**
	 * If the view is updated the new bins are set and then painted.
	 * 
	 * @param bins the new bins to display.
	 */
	public void updateView() {
		setText();
		repaint();
	}
}