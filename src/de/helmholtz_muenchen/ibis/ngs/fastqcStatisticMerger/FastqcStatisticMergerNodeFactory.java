package de.helmholtz_muenchen.ibis.ngs.fastqcStatisticMerger;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FastqcStatisticMergerNodeModel" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FastqcStatisticMergerNodeFactory 
        extends NodeFactory<FastqcStatisticMergerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FastqcStatisticMergerNodeModel createNodeModel() {
        return new FastqcStatisticMergerNodeModel();
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
    public NodeView<FastqcStatisticMergerNodeModel> createNodeView(final int viewIndex,
            final FastqcStatisticMergerNodeModel nodeModel) {
        return new FastqcStatisticMergerNodeView(nodeModel);
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
        return new FastqcStatisticMergerNodeDialog();
    }

}

