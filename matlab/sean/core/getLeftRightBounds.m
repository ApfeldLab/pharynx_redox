function bounds = getLeftRightBounds(segStack)
    bboxes = getBBoxes(segStack);
    bounds = uint8([bboxes(:,1) bboxes(:,1) + bboxes(:, 3)]); % left, right
end