function theta_Fields = convert_theta_to_theta_Fields(theta)

[Fields, GradxFields, GradyFields] = convert_theta_to_Fields(theta);

theta_Fields(:,[1:2],:) = Fields(:,[1:2],:);
theta_Fields(:,[3,5],:) = GradxFields(:,1:2,:);
theta_Fields(:,[4,6],:) = GradyFields(:,1:2,:);
theta_Fields(:,[7:9],:) = Fields(:,[3:5],:);

end
