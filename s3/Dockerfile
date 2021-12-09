FROM alpine

COPY --from=minio/minio /opt/bin/minio /usr/bin/minio
COPY --from=minio/mc /bin/mc /usr/bin/mc

EXPOSE 9000

ADD ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT [ "/entrypoint.sh" ]
