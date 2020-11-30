FROM alpine

COPY --from=minio/minio /usr/bin/minio /usr/bin/minio
COPY --from=minio/mc /usr/bin/mc /usr/bin/mc

EXPOSE 9000

ADD ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT [ "/entrypoint.sh" ]
