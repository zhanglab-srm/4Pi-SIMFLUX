function [phase, amp, offset] = MyCalPhase3_new_parallel(intlist, phase0)
% 并行版本的相位计算函数
% 输入:
%   intlist - N×3 矩阵，每行代表一组 [int1, int2, int3] 数据
%   phase0 - 标量或长度为N的向量，初始相位
% 输出:
%   phase - N×1 向量，计算得到的相位
%   amp - N×1 向量，计算得到的幅度
%   offset - N×1 向量，计算得到的偏移

% 获取数据行数
N = size(intlist, 1);

% 检查 phase0 的维度，如果是标量则扩展为向量
if isscalar(phase0)
    phase0 = repmat(phase0, N, 1);
end

% 提取三列数据
int1 = intlist(:, 1);
int2 = intlist(:, 2);
int3 = intlist(:, 3);

% 向量化计算
cp = cos(phase0);
sp = sin(phase0);

% 计算偏移
offset = (int2 .* cp - (int1 + int3) / 2) ./ (cp - 1);

% 减去偏移
int1 = int1 - offset;
int2 = int2 - offset;
int3 = int3 - offset;

% 计算幅度
amp = sqrt(int2.^2 + (int1 - int2 .* cp).^2 ./ (sp.^2));

% 计算中间变量
sx = int2 ./ amp;
cx = (int1 - int2 .* cp) ./ sp ./ amp;

% 计算临时相位
tp = acos(cx);

% 根据 sx 的符号确定最终相位
phase = zeros(N, 1);
negative_mask = sx < 0;
phase(negative_mask) = -real(tp(negative_mask));
phase(~negative_mask) = real(tp(~negative_mask));

end