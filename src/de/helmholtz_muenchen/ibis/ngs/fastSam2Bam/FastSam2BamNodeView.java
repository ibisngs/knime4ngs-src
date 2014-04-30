package de.helmholtz_muenchen.ibis.ngs.fastSam2Bam;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FastSam2Bam" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FastSam2BamNodeView extends NodeView<FastSam2BamNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FastSam2BamNodeModel})
     */
    protected FastSam2BamNodeView(final FastSam2BamNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() { }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() { }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() { }
}

