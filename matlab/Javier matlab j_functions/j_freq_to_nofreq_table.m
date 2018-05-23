function a = j_freq_to_nofreq_table(strain)

% transforms tables with age of death , freq, censoring
% to age of death, censoring

%strain= {s1, s2};

strain_number = size(strain,2);

for n = 1:strain_number
    strain_nf{n} = [];
    for m = 1:size(strain{n},1)
        for freq = 1:strain{n}(m,2)
            strain_nf{n} = cat(1, strain_nf{n},strain{n}(m,[1 3]));
        end
    end
end

a = strain_nf;




