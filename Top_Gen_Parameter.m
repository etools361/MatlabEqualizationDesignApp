%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 生成设计参数
%--------------------------------------------------------------------------
R = logspace(log10(1.5),log10(10),9);%1:1:9;
NN = 1:10;
mNN = length(NN);
aSol = [];
wSol = [];
deltaSol = [];
xSol     = [];
C = 2:10;
tic;
for jj=1:length(R)
    for ii=1:mNN
        N = NN(ii);
        n = 10^(R(jj)/2);
        m = 10^(-R(jj)/2);
        WC = sqrt(m*n);
        for kk=1:length(C)
            c = C(kk);
            maxIter = 30;
            tol = 1e-13;
            if ii<5
                [aSol{jj,ii,kk},wSol{jj,ii,kk},deltaSol{jj,ii,kk},xSol{jj,ii,kk}] = funRemezEquirippleRational(N, m, n, c, maxIter, tol, 0);
            else
                if mod(N,2)            
                    aSel = [aSol{jj,ii-1,kk}(1);aSol{jj,ii-1,kk}];
                    wSel = logspace(log10(sqrt(n*m)),log10(wSol{jj,ii-1,kk}(end)),floor((N+1)/2))';
                else
                    aSel = [aSol{jj,ii-1,kk}(1);aSol{jj,ii-1,kk}(2:end)];
                    wSel = [sqrt(sqrt(n*m)*wSol{jj,ii-1,kk}(2));wSol{jj,ii-1,kk}(2:end)];
                end
                [aSol{jj,ii,kk},wSol{jj,ii,kk},deltaSol{jj,ii,kk},xSol{jj,ii,kk}] = funRemezEquirippleRational(N, m, n, c, maxIter, tol, 0, abs(aSel), abs(wSel));
            end
        end
    end
end
toc;
wSol2 = [];
for ii=1:length(NN)
    for jj=1:length(R)
        for kk=1:length(C)
        wSol2{jj,ii,kk} = log10(wSol{jj,ii,kk})./(R(jj)/2);
%         plot(wSol2{jj,ii,kk}, aSol{jj,ii,kk}, '-*', 'linewidth', 2);
%         grid on;
%         hold on;
        end
    end
end
% hold off;

wSeg = [];
aSeg = [];
for ii=1:length(NN)
    for jj=1:length(R)
        for kk=1:length(wSol2{jj,ii,kk})
            for mm=1:length(C)
                wSeg{ii}(jj,kk,mm) = wSol2{jj,ii,mm}(kk);
                aSeg{ii}(jj,kk,mm) = aSol{jj,ii,mm}(kk);
            end
        end
    end
end

% for ii=8
%     plot(R, wSeg{ii}, '-*', 'linewidth', 2);
%     hold on;
%     grid on;
% end
% hold off;

fprintf('wTabNxCx = {\n');
for mm=1:length(C)
    fprintf('{\n')
    for ii=1:length(NN)
        [a, b] = size(wSeg{ii}(:,:,mm));
    %     fprintf('-------%d---------\n', R(ii));
    %     fprintf('wTabN%dC%d = [\n', NN(ii), C(1));
        fprintf('[\n');
        for jj=1:b
            for kk=1:a
                if kk==a
                    fprintf('%0.3f;\n', wSeg{ii}(kk,jj,mm));
                else
                    fprintf('%0.3f,', wSeg{ii}(kk,jj,mm));
                end
            end
        end
        fprintf('];\n');
        fprintf('[\n');
        for jj=1:b
            for kk=1:a
                if kk==a
                    fprintf('%0.3f;\n', aSeg{ii}(kk,jj,mm));
                else
                    fprintf('%0.3f,', aSeg{ii}(kk,jj,mm));
                end
            end
        end
        fprintf('];\n');
    end
    fprintf('};\n');
end
fprintf('};\n');



% [ax2,wx2,deltaSolx2,xSolx2] = funRemezEquirippleRational(Nx, m1, n1, c1, maxIter, tol, 7, axc0', 10.^(wxc0'*Rx/2));



