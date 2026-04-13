get_num_of_dims(joint::LinearJoint) = get_num_of_dims(joint.body)
get_num_of_dims(joint::FixedPointJoint) = get_num_of_dims(joint.body)
