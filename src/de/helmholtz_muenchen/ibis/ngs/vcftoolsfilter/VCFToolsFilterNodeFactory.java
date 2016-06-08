package de.helmholtz_muenchen.ibis.ngs.vcftoolsfilter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "LOFFilter" Node.
 * 
 *
 * @author tim.jeske
 */
public class VCFToolsFilterNodeFactory 
        extends NodeFactory<VCFToolsFilterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VCFToolsFilterNodeModel createNodeModel() {
        return new VCFToolsFilterNodeModel();
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
    public NodeView<VCFToolsFilterNodeModel> createNodeView(final int viewIndex,
            final VCFToolsFilterNodeModel nodeModel) {
        return new HTENodeView<VCFToolsFilterNodeModel>(nodeModel);
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
        return new VCFToolsFilterNodeDialog();
    }

}

