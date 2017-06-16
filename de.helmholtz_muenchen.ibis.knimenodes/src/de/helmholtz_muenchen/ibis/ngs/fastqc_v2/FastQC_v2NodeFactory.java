package de.helmholtz_muenchen.ibis.ngs.fastqc_v2;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.ngs.fastqc.FastQCNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "FastQC_v2" Node.
 * 
 *
 * @author Paul Hager
 */
public class FastQC_v2NodeFactory 
        extends NodeFactory<FastQC_v2NodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FastQC_v2NodeModel createNodeModel() {
        return new FastQC_v2NodeModel();
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
    public NodeView<FastQC_v2NodeModel> createNodeView(final int viewIndex,
            final FastQC_v2NodeModel nodeModel) {
        return new HTENodeView<FastQC_v2NodeModel>(nodeModel);
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
        return new FastQC_v2NodeDialog();
    }

}

