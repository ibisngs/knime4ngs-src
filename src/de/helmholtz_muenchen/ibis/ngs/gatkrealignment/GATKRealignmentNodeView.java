package de.helmholtz_muenchen.ibis.ngs.gatkrealignment;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "GATKRealignment" Node.
 * 
 *
 * @author 
 */
public class GATKRealignmentNodeView extends NodeView<GATKRealignmentNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link GATKRealignmentNodeModel})
     */
    protected GATKRealignmentNodeView(final GATKRealignmentNodeModel nodeModel) {
        super(nodeModel);

        //instantiate the components of the view here.

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {

        // retrieve the new model from your nodemodel and 
        // update the view.
        GATKRealignmentNodeModel nodeModel = 
            (GATKRealignmentNodeModel)getNodeModel();
        assert nodeModel != null;
        
        // be aware of a possibly not executed nodeModel! The data you retrieve
        // from your nodemodel could be null, emtpy, or invalid in any kind.
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    
        // things to do when closing the view
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {

        // things to do when opening the view
    }

}

