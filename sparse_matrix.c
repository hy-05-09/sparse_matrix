
#include "sparseMTran.h"

extern int allocSize;      // size of allocated memory currently used

matrixPtr mread(void) {
    matrixPtr node = NULL;      
    matrixPtr* hdnode = NULL;   
    int nheads = 0;
    int nrows, ncols, nterms, i;
    int row, col, val, currow;
    matrixPtr tmp, last;

    scanf("%d %d %d", &nrows, &ncols, &nterms);
    nheads = (ncols > nrows) ? ncols : nrows; //row와 col 중에 더 큰 수 할당
    node = newNode();
    allocSize += sizeof(matrixNode); //node의 메모리 크기
    node->tag = entry;
    node->u.entry.row = nrows;
    node->u.entry.col = ncols;
    node->u.entry.val = nterms;

    if (!nheads) { //nHeads=0일 경우
        node->right = node; //linked list 생성 x
        return node;
    }

    MALLOC(hdnode, matrixPtr*, nheads * sizeof(matrixPtr)); //헤더노드 포인터 malloc
    allocSize += nheads * sizeof(matrixPtr); //헤더노드 포인터 메모리 크기
    for (i = 0; i < nheads; i++) {
        tmp = newNode();
        allocSize += sizeof(matrixNode); //header node 메모리 크기
        hdnode[i] = tmp;
        hdnode[i]->tag = head;
        hdnode[i]->right = tmp;
        hdnode[i]->u.next = tmp;
    }
    currow = 0; //현재 row 
    last = hdnode[0];

    for (i = 0; i < nterms; i++) {
        scanf("%d %d %d", &row, &col, &val);
        if (row > currow) { //tmp의 row가 현재 row보다 클 때
            last->right = hdnode[currow]; //이전 tmp의 right에 헤더노드 연결
            currow = row; //현재 row 갱신
            last = hdnode[row]; //last 갱신
        }
        tmp = newNode();
        allocSize += sizeof(matrixNode); //element node 메모리 크기
        tmp->u.entry.row = row;
        tmp->u.entry.col = col;
        tmp->u.entry.val = val;
        last->right = tmp; //이전 tmp의 right에 현재 tmp 연결
        last = tmp; //last 갱신
        hdnode[col]->u.next->down = tmp; //해당 행의 가장 아래였던 node의 down에 tmp 연결
        hdnode[col]->u.next = tmp; //해당 행 가장 아래 node 갱신
    }

    last->right = hdnode[currow]; //해당 열 마지막 node에 헤더노드 연결
    for (i = 0; i < ncols; i++)
        hdnode[i]->u.next->down = hdnode[i]; //해당 행 마지막 node에 헤더노드 연결
    for (i = 0; i < nheads - 1; i++)
        hdnode[i]->u.next = hdnode[i + 1]; //헤더노드 next에 다음 헤더노드 연결
    hdnode[nheads - 1]->u.next = node; //마지막 헤더노드에 node 연결
    node->right = hdnode[0]; //node의 right에 첫번째 헤더노드 연결

   						 
    allocSize -= nheads * sizeof(matrixPtr); //헤더노드 포인터 메모리 해제
    free(hdnode);	// Deallocate hdnode array(not needed)

    return node;
}

void merase(matrixPtr* node) {

    matrixPtr x, y, head = (*node)->right;
    for (int i = 0; i < (*node)->u.entry.row; i++) {
        y = head->right;
        while (y != head) { //해당 열에 node가 남아있을 때
            x = y; y = y->right; 
            allocSize -= sizeof(matrixNode); //노드 메모리 크기 빼기
            free(x); //메모리 해제
        }
        x = head; head = head->u.next; allocSize -= sizeof(matrixNode); free(x); //헤더노드 해제 및 메모리 계산
    }
    y = head;
    while (y != *node) { //col이 row보다 클 때
        x = y; y = y->u.next; allocSize -= sizeof(matrixNode); free(x); //남은 메모리 해제
    }
    allocSize -= sizeof(matrixNode); //node 메모리 계산
    free(*node);

	*node = NULL;
	
    return;
}

matrixPtr mtranspose(matrixPtr A) {
    
    matrixPtr node = NULL;  // header pointer for the transpose of A
    matrixPtr* hdnode = NULL;   // head node array for the transpose of A
    int nHeads = MAX2(A->u.entry.col, A->u.entry.row);  // # of head nodes
    int nrows, ncols, nterms, i;
    int row, col, val, curcol;
    matrixPtr tmp, last;
    matrixPtr x,y,head2 = A->right;

    nrows = A->u.entry.col; //A의 열 -> 새로운 노드의 행으로 저장
    ncols = A->u.entry.row; //A의 행 -> 새로운 노드의 열로 저장
    nterms = A->u.entry.val;

    node = newNode();
    allocSize += sizeof(matrixNode); //node 메모리 계산
    node->tag = entry;
    node->u.entry.row = nrows;
    node->u.entry.col = ncols;
    node->u.entry.val = nterms;

    if (!nHeads) { //nHeads==0일 때
        node->right = node;
        return node;
    }

    MALLOC(hdnode, matrixPtr*, nHeads * sizeof(matrixPtr)); //헤더노드 포인터 할당
    allocSize += nHeads * sizeof(matrixPtr); //헤더노드 포인터 메모리
    for (i = 0; i < nHeads; i++) {
        tmp = newNode();
        allocSize += sizeof(matrixNode); //헤더노드 메모리
        hdnode[i] = tmp;
        hdnode[i]->tag = head;
        hdnode[i]->right = tmp;
        hdnode[i]->u.next = tmp;
    }
    curcol = 0; //A의 현재 col
    last = hdnode[0];

    for (i = 0; i < A->u.entry.col; i++) { //A는 col 기준으로
        y = head2->down; //col 기준이므로 down
        while (y != head2) { //해당 col에 node가 남아있을 때
            row = y->u.entry.col; //바꿔 저장
            col = y->u.entry.row; //바꿔 저장
            val = y->u.entry.val;

            if (row > curcol) { //row(A의 col)이 curcol보다 클 때
                last->right = hdnode[curcol]; //새로운 노드의 last의 right에 헤더노드 연결
                curcol = row; //A의 curcol 갱신
                last = hdnode[row]; //새로운 노드의 last 갱신
            }
            tmp = newNode();
            allocSize += sizeof(matrixNode); //새로운 노드 메모리
            tmp->tag = entry;
            tmp->u.entry.row = row;
            tmp->u.entry.col = col;
            tmp->u.entry.val = val;
            last->right = tmp; 
            last = tmp;
            hdnode[col]->u.next->down = tmp; //새로운 노드의 행 포인터 연결
            hdnode[col]->u.next = tmp; //새로운 노드의 행 포인터 갱신
            y = y->down;   //A 다음 노드로 이동
        }
        head2 = head2->u.next; //A 다음 행으로 이동
    }
   

    last->right = hdnode[curcol];
    for (i = 0; i < ncols; i++)
        hdnode[i]->u.next->down = hdnode[i];
    for (i = 0; i < nHeads - 1; i++)
        hdnode[i]->u.next = hdnode[i + 1];
    hdnode[nHeads - 1]->u.next = node;
    node->right = hdnode[0];

   
    allocSize -= nHeads * sizeof(matrixPtr); //헤더노드 포인터 메모리
    free(hdnode); //헤더노드 포인터 해제

    return node;
}

void mwriteFull(matrixPtr node) {
    int k, m, j=0;
    matrixPtr tmp, head = node->right;
    printf("[%d X %d]\n", node->u.entry.row, node->u.entry.col);
    tmp = head->right;
    for (int i = 0; i < node->u.entry.row; i++) { //열만큼 실행
        k = i, m = 0;
      

        while (true) {
            if (m >= node->u.entry.col) { //행 수가 다 찼을 때
                break; //다음 열로
            }
            if (tmp->u.entry.row == k) { //행이 같을 때
                if (tmp->u.entry.col == m) { //열이 같을 때
                    if (m<(node->u.entry.col-1)) //마지막 행이 아닐 때
                        printf("%5d ", tmp->u.entry.val); //출력
                    else if (m==(node->u.entry.col-1)) //마지막 행일 때 
                        printf("%5d\n", tmp->u.entry.val); //출력
                    m++; //행++
                    j++; //출력할 val++
                    tmp = tmp->right; //다음 노드로 이동
                    if (tmp == head){  //해당 열의 노드 끝
                        while (m < node->u.entry.col) { //행 수 다 안 찼을 때
                            if (m < (node->u.entry.col - 1))
                                printf("%5d ", 0); //0 출력
                            else if (m == (node->u.entry.col - 1))
                                printf("%5d\n", 0); //0 출력
                            m++; //행++
                        }
                        head = head->u.next; //헤더 이동
                        tmp = head->right; //노드 이동
                        break; //다음 열 실행
                    }
                    continue; //다음 행 실행
                }
                else if (tmp->u.entry.col > m) { //현재 행이 노드의 행보다 작을 때
                    if (m < (node->u.entry.col - 1))
                        printf("%5d ", 0); //0 출력
                    else if (m == (node->u.entry.col - 1))
                        printf("%5d\n", 0); //0 출력
                    m++; //행++
                    continue;
                }
            }
            else if (tmp==head) { //해당 열 노드 출력 끝
                while (m != node->u.entry.col) { //행 안 끝났을 때
                    if (m < (node->u.entry.col - 1))
                        printf("%5d ", 0); //0 출력
                    else if (m == (node->u.entry.col - 1))
                        printf("%5d\n", 0); //0 출력
                    m++; //행++
                }
                head = head->u.next; //헤더 이동
                tmp = head; 
                tmp = tmp->right; //노드 이동
                break;
            }
        }
    }
    
    return;
}