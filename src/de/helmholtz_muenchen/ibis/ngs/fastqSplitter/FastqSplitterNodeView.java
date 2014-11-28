package de.helmholtz_muenchen.ibis.ngs.fastqSplitter;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FastqSplitter" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FastqSplitterNodeView extends NodeView<FastqSplitterNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FastqSplitterNodeModel})
     */
    protected FastqSplitterNodeView(final FastqSplitterNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
    }

}

