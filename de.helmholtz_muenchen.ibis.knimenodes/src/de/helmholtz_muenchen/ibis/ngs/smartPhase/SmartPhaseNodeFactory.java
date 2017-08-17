package de.helmholtz_muenchen.ibis.ngs.smartPhase;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "SmartPhase" Node.
 * KNIME integration of the smartPhase program that phases filtered rare variants in desired genomic regions to determine possible compound heterozygosity.
 *
 * @author Paul Hager
 */
public class SmartPhaseNodeFactory 
        extends NodeFactory<SmartPhaseNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SmartPhaseNodeModel createNodeModel() {
        return new SmartPhaseNodeModel();
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
    public NodeView<SmartPhaseNodeModel> createNodeView(final int viewIndex,
            final SmartPhaseNodeModel nodeModel) {
        return new HTENodeView<SmartPhaseNodeModel>(nodeModel);
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
        return new SmartPhaseNodeDialog();
    }

}

