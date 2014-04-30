package de.helmholtz_muenchen.ibis.ngs.fastSam2Bam;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FastSam2Bam" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FastSam2BamNodeFactory 
        extends NodeFactory<FastSam2BamNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FastSam2BamNodeModel createNodeModel() {
        return new FastSam2BamNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<FastSam2BamNodeModel> createNodeView(final int viewIndex,
            final FastSam2BamNodeModel nodeModel) {
        return new FastSam2BamNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new FastSam2BamNodeDialog();
    }

}

